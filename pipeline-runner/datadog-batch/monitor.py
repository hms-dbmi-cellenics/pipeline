from datetime import datetime
import os
from time import sleep
from datadog_api_client import ApiClient, Configuration
from datadog_api_client.v2.api.metrics_api import MetricsApi
from datadog_api_client.v2.model.metric_payload import MetricPayload
from datadog_api_client.v2.model.metric_point import MetricPoint
from datadog_api_client.v2.model.metric_resource import MetricResource
from datadog_api_client.v2.model.metric_series import MetricSeries

from instance_metadata import get_instance_metadata
from cgroup_v1_stats import CGroupV1Stats
from cgroup_v2_stats import CGroupV2Stats

from constants import DatadogMetricTypes


METRICS_IN_SERIES = 15
COLLECTION_INTERVAL = 1 #s

def is_cgroup_v2():
    # Check as suggested by: https://unix.stackexchange.com/a/668244/553416
    return os.path.exists("/sys/fs/cgroup/cgroup.controllers")

def collect_metrics(is_cgroup_v2):

    stats = {}
    if is_cgroup_v2:
        collector = CGroupV2Stats()
        stats = collector.collect()
    else:
        collector = CGroupV1Stats()
        stats = collector.collect()

    return {
        "timestamp": int(datetime.now().timestamp()),
        **stats
    }


def format_metrics(metrics, resources, tags):

    payload = MetricPayload(
        series=[
            MetricSeries(
                metric=body['metric'],
                type=body['type'],
                unit=body['unit'],
                points=[
                    MetricPoint(
                        timestamp=metric["timestamp"],
                        value=metric[key],
                    ) for metric in metrics
                ],
                resources=resources,
                tags=tags
            ) for key, body in DatadogMetricTypes.items()
        ]
    )

    return payload

if __name__ == "__main__":

    print("Datadog subprocess launched")

    print('*** os.getenv("DD_API_KEY")', os.getenv("DD_API_KEY"))
    print('*** os.getenv("DD_APP_KEY")', os.getenv("DD_APP_KEY"))
    print('*** os.getenv("AWS_DEFAULT_REGION")', os.getenv("AWS_DEFAULT_REGION"))
    print('*** os.getenv("ACTIVITY_ARN")', os.getenv("ACTIVITY_ARN"))
    print('*** os.getenv("CLUSTER_ENV")', os.getenv("CLUSTER_ENV"))
    print('*** os.getenv("AWS_ACCOUNT_ID")', os.getenv("AWS_ACCOUNT_ID"))
    print('*** os.getenv("SANDBOX_ID")', os.getenv("SANDBOX_ID"))

    configuration = Configuration()
    configuration.api_key["apiKeyAuth"] = os.getenv("DD_API_KEY")
    configuration.api_key["appKeyAuth"] = os.getenv("DD_APP_KEY")
    configuration.server_variables["site"] = "datadoghq.eu"

    instance_meta = get_instance_metadata([
        "instance-id"
        "instance-type"
        "hostname"
    ])

    print("*** instance meta", instance_meta)
    print("*** configuration.api_key['apiKeyAuth']", configuration.api_key["apiKeyAuth"])
    print("*** configuration.api_key['appKeyAuth']", configuration.api_key["appKeyAuth"])

    resources = [
        MetricResource(
            name=f"{instance_meta.get('instance-id', 'batch-instance')}",
            type="host",
        )
    ]

    tags = [
        f"instanceId:{instance_meta.get('instance-id', 'batch-instance')}",
        f"instanceType:{instance_meta.get('instance-type', '')}",
        f"hostname:{instance_meta.get('hostname', '')}",
        f"service:batch",
        f"region:{os.getenv('AWS_DEFAULT_REGION')}",
        f"activityId:{os.getenv('ACTIVITY_ARN')}",
        f"env:{os.getenv('CLUSTER_ENV')}",
        f"awsAccountId:{os.getenv('AWS_ACCOUNT_ID')}",
        f"sandboxId:{os.getenv('SANDBOX_ID')}",
    ]

    print("*** is_cgroup_v2", is_cgroup_v2())
    print("*** tags", tags)

    print(f"Reporting metrics to Datadog every {METRICS_IN_SERIES}s")

    with ApiClient(configuration) as api_client:
        api_instance = MetricsApi(api_client)

        while True:
            series = []

            i = 0
            while i < METRICS_IN_SERIES:
                series.append(collect_metrics(is_cgroup_v2()))
                i += 1
                sleep(COLLECTION_INTERVAL)

            payload = format_metrics(series, resources, tags)

            print("*** payload", payload)

            response = api_instance.submit_metrics(body=payload)

            print("*** response", response)

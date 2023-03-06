from datetime import datetime
import os
from math import floor
from time import sleep
from datadog_api_client import ApiClient, Configuration
from datadog_api_client.v2.api.metrics_api import MetricsApi
from datadog_api_client.v2.model.metric_resource import MetricResource

from instance_metadata import get_instance_metadata
from cgroup_v1_stats import CGroupV1Stats
from cgroup_v2_stats import CGroupV2Stats

from format_metrics import format_metrics


REPORTING_INTERVAL = 15 #s
COLLECTION_INTERVAL = 1 #s
NUM_METRICS_IN_SERIES = floor(REPORTING_INTERVAL / COLLECTION_INTERVAL)

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


if __name__ == "__main__":
    print("Datadog monitoring subprocess launched")

    configuration = Configuration()
    configuration.api_key["apiKeyAuth"] = os.getenv("DD_API_KEY")
    configuration.api_key["appKeyAuth"] = os.getenv("DD_APP_KEY")
    configuration.server_variables["site"] = "datadoghq.eu"

    instance_meta = get_instance_metadata([
        "instance-id",
        "instance-type",
        "hostname"
    ])

    resources = [
        MetricResource(
            name=f"{instance_meta.get('instance-id', '')}",
            type="host",
        )
    ]

    tags = [
        f"instanceId:{instance_meta.get('instance-id', '')}",
        f"instanceType:{instance_meta.get('instance-type', '')}",
        f"hostname:{instance_meta.get('hostname', '')}",
        f"service:batch",
        f"region:{os.getenv('AWS_DEFAULT_REGION')}",
        f"experimentId:{os.getenv('EXPERIMENT_ID')}",
        f"activityId:{os.getenv('ACTIVITY_ARN')}",
        f"env:{os.getenv('CLUSTER_ENV')}",
        f"awsAccountId:{os.getenv('AWS_ACCOUNT_ID')}",
        f"sandboxId:{os.getenv('SANDBOX_ID')}",
    ]

    print(f"Collecting data every {COLLECTION_INTERVAL}s, reporting to Datadog every {REPORTING_INTERVAL}s")

    with ApiClient(configuration) as api_client:
        api_instance = MetricsApi(api_client)

        while True:
            series = []

            i = 0
            while i < NUM_METRICS_IN_SERIES:
                series.append(collect_metrics(is_cgroup_v2()))
                i += 1
                sleep(COLLECTION_INTERVAL)

            payload = format_metrics(series, resources, tags)
            response = api_instance.submit_metrics(body=payload)

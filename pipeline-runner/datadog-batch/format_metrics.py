from datadog_api_client.v2.model.metric_payload import MetricPayload
from datadog_api_client.v2.model.metric_point import MetricPoint
from datadog_api_client.v2.model.metric_series import MetricSeries

from constants import ReportedDatadogMetrics, MetricsGroup

def get_cpu_metric_series(series, resources, tags):
    metric_series = {}
    metric_group = MetricsGroup.CPU
    reported_metrics = ReportedDatadogMetrics[MetricsGroup.CPU]

    for time_step in series:
        for metric_name in reported_metrics.keys():
            if not metric_series.get(metric_name):
                metric_series[metric_name] = []

            metric_series[metric_name].append(
                MetricPoint(
                    timestamp=time_step["timestamp"],
                    value=time_step[metric_group][metric_name],
                )
            )

    result = []
    for metric_name, body in reported_metrics.items():
        result.append(
            MetricSeries(
                metric=body['metric'],
                type=body['type'],
                unit=body['unit'],
                points=metric_series[metric_name],
                resources=resources,
                tags=tags
            )
        )

    return result


def get_memory_metric_series(series, resources, tags):
    metric_series = {}
    metric_group = MetricsGroup.MEMORY
    reported_metrics = ReportedDatadogMetrics[metric_group]

    for time_step in series:
        for metric_name in reported_metrics.keys():
            if not metric_series.get(metric_name):
                metric_series[metric_name] = []

            metric_series[metric_name].append(
                MetricPoint(
                    timestamp=time_step["timestamp"],
                    value=time_step[metric_group][metric_name],
                )
            )

    result = []
    for metric_name, body in reported_metrics.items():
        result.append(
            MetricSeries(
                metric=body['metric'],
                type=body['type'],
                unit=body['unit'],
                points=metric_series[metric_name],
                resources=resources,
                tags=tags
            )
        )

    return result


def get_io_metric_series(series, resources, tags):
    metric_series = {}
    metric_group = MetricsGroup.IO
    reported_metrics = ReportedDatadogMetrics[metric_group]

    for time_step in series:
        for device_id, device_metric in time_step[metric_group].items():
            for metric_name in reported_metrics.keys():
                if not metric_series.get(device_id):
                    metric_series[device_id] = {}

                if not metric_series[device_id].get(metric_name):
                    metric_series[device_id][metric_name] = []

                metric_series[device_id][metric_name].append(
                    MetricPoint(
                        timestamp=time_step["timestamp"],
                        value=device_metric[metric_name],
                    )
                )

    result = []
    for metric_name, body in reported_metrics.items():
        for device_id, device_metric in metric_series.items():
            result.append(
                MetricSeries(
                    metric=body['metric'],
                    type=body['type'],
                    unit=body['unit'],
                    points=device_metric[metric_name],
                    resources=resources,
                    tags=[*tags, f"deviceId:{device_id}"]
                )
            )

    return result


def format_metrics(series, resources, tags):

    cpu_metric_series = get_cpu_metric_series(series, resources, tags)
    memory_metric_series = get_memory_metric_series(series, resources, tags)
    io_metric_series = get_io_metric_series(series, resources, tags)

    return MetricPayload(series=[
        *cpu_metric_series,
        *memory_metric_series,
        *io_metric_series
    ])
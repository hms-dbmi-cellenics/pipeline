from enum import Enum

from datadog_api_client.v2.model.metric_intake_type import MetricIntakeType

DatadogMetrics = Enum(
    'DatadogMetrics',
    [
        'CONTAINER_CPU_USAGE',
        'CONTAINER_CPU_USER',
        'CONTAINER_CPU_SYSTEM',
        'CONTAINER_CPU_THROTTLED',
        'CONTAINER_CPU_THROTTLED_PERIODS',
        'CONTAINER_MEMORY_USAGE',
        'CONTAINER_MEMORY_LIMIT',
        'CONTAINER_MEMORY_SWAP',
        'CONTAINER_MEMORY_CACHE',
        'CONTAINER_MEMORY_RSS',
        'CONTAINER_MEMORY_KERNEL',
        'CONTAINER_MEMORY_OOM_EVENTS',
        'CONTAINER_IO_READ',
        'CONTAINER_IO_WRITE',
        'CONTAINER_IO_READ_OPERATIONS',
        'CONTAINER_IO_WRITE_OPERATIONS',
    ]
)

MetricsGroup = Enum(
    'DatadogMetricsGroup',
    [
        "CPU",
        "MEMORY",
        "IO"
    ]
)

# For metric name, type and unit, refer to https://datadoghq.eu/metric/summary
ReportedDatadogMetrics = {
    MetricsGroup.CPU: {
        DatadogMetrics.CONTAINER_CPU_USAGE: {
            "metric": "container.cpu.usage",
            "type": MetricIntakeType.RATE,
            "unit": "nanosecond"
        },
        DatadogMetrics.CONTAINER_CPU_USER: {
            "metric": "container.cpu.user",
            "type": MetricIntakeType.RATE,
            "unit": "nanosecond"
        },
        DatadogMetrics.CONTAINER_CPU_SYSTEM: {
            "metric": "container.cpu.system",
            "type": MetricIntakeType.RATE,
            "unit": "nanosecond"
        },
        DatadogMetrics.CONTAINER_CPU_THROTTLED: {
            "metric": "container.cpu.throttled",
            "type": MetricIntakeType.RATE,
            "unit": "nanosecond"
        },
        DatadogMetrics.CONTAINER_CPU_THROTTLED_PERIODS: {
            "metric": "container.cpu.throttled.periods",
            "type": MetricIntakeType.RATE,
            "unit": "nanosecond"
        },
    },
    MetricsGroup.MEMORY: {
        DatadogMetrics.CONTAINER_MEMORY_USAGE: {
            "metric": "container.memory.usage",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_LIMIT: {
            "metric": "container.memory.limit",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_SWAP: {
            "metric": "container.memory.swap",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_CACHE: {
            "metric": "container.memory.cache",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_RSS: {
            "metric": "container.memory.rss",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_KERNEL: {
            "metric": "container.memory.kernel",
            "type": MetricIntakeType.GAUGE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_MEMORY_OOM_EVENTS: {
            "metric": "container.memory.oom_events",
            "type": MetricIntakeType.GAUGE,
            "unit": None
        },
    },
    MetricsGroup.IO: {
        DatadogMetrics.CONTAINER_IO_READ: {
            "metric": "container.io.read",
            "type": MetricIntakeType.RATE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_IO_WRITE: {
            "metric": "container.io.write",
            "type": MetricIntakeType.RATE,
            "unit": "byte"
        },
        DatadogMetrics.CONTAINER_IO_READ_OPERATIONS: {
            "metric": "container.io.read.operations",
            "type": MetricIntakeType.RATE,
            "unit": None
        },
        DatadogMetrics.CONTAINER_IO_WRITE_OPERATIONS: {
            "metric": "container.io.write.operations",
            "type": MetricIntakeType.RATE,
            "unit": None
        },
    }
}
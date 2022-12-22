# These are compiled from:
# https://github.com/DataDog/datadog-agent/blob/main/pkg/collector/corechecks/containers/generic/processor.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv1_cpu.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv1_memory.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv1_io.go

from constants import DatadogMetrics, MetricsGroup

BASE_PATH = '/sys/fs/cgroup'
CPU_PATH = f'{BASE_PATH}/cpu'
MEMORY_PATH = f'{BASE_PATH}/memory'
IO_PATH = f'{BASE_PATH}/blkio'

class CGroupV1Stats:

    def __init__(self):
        self.stats = {
            MetricsGroup.CPU: {},
            MetricsGroup.MEMORY: {},
            MetricsGroup.IO: {}
        }
        self._read_cpu_stats()
        self._read_memory_stats()
        self._read_io_stats()

    def _read_float(self, folder, file):
        with open(f"{folder}/{file}") as f:
            return float(f.readline())

    def _read_multi_stats(self, folder, file, group, mapping):
        lines = []
        value_mapping = {}

        with open(f"{folder}/{file}") as f:
            lines = f.readlines()

        for line in lines:
            key, value = line.split(" ")
            value_mapping[key] = value

        for stats_key, value_key in mapping.items():
            self.stats[stats_key] = float(value_mapping[value_key])

    def _parse_io_stats(self, folder, file, group, io_read_key, io_write_key):
        # Example of IO file
        # 259:0 Read 38
        # 259:0 Write 2
        # 259:0 Sync 40
        # 259:0 Async 0
        # 259:0 Discard 0
        # 259:0 Total 40
        # Total 40

        with open(f"{folder}/{file}") as f:
            lines = f.readlines()

            for line in lines:
                # Lines with less than 3 fields do not contain device-level information
                if len(line) < 3:
                    continue

                device_id, key, value = line.split(" ")

                if not self.stats[group].get(device_id):
                    self.stats[group][device_id] = {}

                if key == "Read":
                    self.stats[group][device_id][io_read_key] = float(value)
                if key == "Write":
                    self.stats[group][device_id][io_write_key] = float(value)

    def _read_cpu_stats(self):
        self.stats[MetricsGroup.CPU][DatadogMetrics.CONTAINER_CPU_USAGE] = self._read_float(CPU_PATH, "cpuacct.usage")

        self._read_multi_stats(CPU_PATH, "cpuacct.stat", MetricsGroup.CPU,
        {
            DatadogMetrics.CONTAINER_CPU_USER: "user",
            DatadogMetrics.CONTAINER_CPU_SYSTEM: "system"
        })

        self._read_multi_stats(CPU_PATH, "cpu.stat", MetricsGroup.CPU,
        {
            DatadogMetrics.CONTAINER_CPU_THROTTLED_PERIODS : "nr_throttled",
            DatadogMetrics.CONTAINER_CPU_THROTTLED: "throttled_time"
        })

    def _read_memory_stats(self):
        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_USAGE] = self._read_float(MEMORY_PATH, 'memory.usage_in_bytes')
        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_KERNEL] = self._read_float(MEMORY_PATH, 'memory.kmem.usage_in_bytes')

        self._read_multi_stats(MEMORY_PATH, 'memory.stat', MetricsGroup.MEMORY,
        {
            DatadogMetrics.CONTAINER_MEMORY_LIMIT : "hierarchical_memory_limit",
            DatadogMetrics.CONTAINER_MEMORY_RSS : "rss",
            DatadogMetrics.CONTAINER_MEMORY_CACHE : "total_cache",
            DatadogMetrics.CONTAINER_MEMORY_SWAP : "total_swap",
        })

        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_OOM_EVENTS] = float(self._read_float(MEMORY_PATH, 'memory.failcnt'))

    def _read_io_stats(self):
        self.stats[MetricsGroup.IO] = {}
        self._parse_io_stats(
            IO_PATH,
            "blkio.throttle.io_service_bytes",
            MetricsGroup.IO,
            DatadogMetrics.CONTAINER_IO_READ,
            DatadogMetrics.CONTAINER_IO_WRITE
        )

        self._parse_io_stats(
            IO_PATH,
            "blkio.throttle.io_serviced",
            MetricsGroup.IO,
            DatadogMetrics.CONTAINER_IO_READ_OPERATIONS,
            DatadogMetrics.CONTAINER_IO_WRITE_OPERATIONS
        )

    def collect(self):
        return self.stats
# These are compiled from:
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv1_cpu.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv1_memory.go

from constants import DatadogMetrics

BASE_PATH = '/sys/fs/cgroup'
CPU_PATH = f'{BASE_PATH}/cpu'
MEMORY_PATH = f'{BASE_PATH}/memory'

class CGroupV1Stats:

    def __init__(self):
        self.stats = {}
        self._read_cpu_stats()
        self._read_memory_stats()

    def _read_float(self, folder, file):
        with open(f"{folder}/{file}") as f:
            return float(f.readline())

    def _read_multi_stats(self, folder, file, mapping):
        lines = []
        value_mapping = {}

        with open(f"{folder}/{file}") as f:
            lines = f.readlines()

        for line in lines:
            key, value = line.split(" ")
            value_mapping[key] = value

        for stats_key, value_key in mapping.items():
            self.stats[stats_key] = float(value_mapping[value_key])

    def _read_cpu_stats(self):
        self.stats[DatadogMetrics.CONTAINER_CPU_USAGE] = self._read_float(CPU_PATH, "cpuacct.usage")

        self._read_multi_stats(CPU_PATH, "cpuacct.stat",
        {
            DatadogMetrics.CONTAINER_CPU_USER: "user",
            DatadogMetrics.CONTAINER_CPU_SYSTEM: "system"
        })

        self._read_multi_stats(CPU_PATH, "cpu.stat", {
            DatadogMetrics.CONTAINER_CPU_THROTTLED_PERIODS : "nr_throttled",
            DatadogMetrics.CONTAINER_CPU_THROTTLED: "throttled_time"
        })

    def _read_memory_stats(self):
        self.stats[DatadogMetrics.CONTAINER_MEMORY_USAGE] = self._read_float(MEMORY_PATH, 'memory.usage_in_bytes')
        self.stats[DatadogMetrics.CONTAINER_MEMORY_KERNEL] = self._read_float(MEMORY_PATH, 'memory.kmem.usage_in_bytes')

        self._read_multi_stats(MEMORY_PATH, 'memory.stat', {
            DatadogMetrics.CONTAINER_MEMORY_LIMIT : "hierarchical_memory_limit",
            DatadogMetrics.CONTAINER_MEMORY_RSS : "rss",
            DatadogMetrics.CONTAINER_MEMORY_CACHE : "total_cache",
            DatadogMetrics.CONTAINER_MEMORY_SWAP : "total_swap",
        })

        self.stats[DatadogMetrics.CONTAINER_MEMORY_OOM_EVENTS] = float(self._read_float(MEMORY_PATH, 'memory.failcnt'))

    def collect(self):
        return self.stats
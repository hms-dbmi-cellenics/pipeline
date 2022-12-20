# These are compiled from:
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv2_cpu.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv2_memory.go

from constants import DatadogMetrics

BASE_PATH = '/sys/fs/cgroup'
MICRO_TO_NANOSECONDS = 1000

class CGroupV2Stats:

    def __init__(self):
        self.stats = {}
        self._read_cpu_stats()
        self._read_memory_stats()

    def _read(self, folder, file):
        with open(f"{folder}/{file}") as f:
            return f.readline()

    def _read_float(self, folder, file):
        return float(self._read(folder, file))

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
        self._read_multi_stats(BASE_PATH, "cpu.stat", {
            DatadogMetrics.CONTAINER_CPU_USAGE: 'usage_usec',
            DatadogMetrics.CONTAINER_CPU_USER: 'user_usec',
            DatadogMetrics.CONTAINER_CPU_SYSTEM: 'system_usec',
            DatadogMetrics.CONTAINER_CPU_THROTTLED_PERIODS: 'nr_throttled',
            DatadogMetrics.CONTAINER_CPU_THROTTLED: 'throttled_usec'
        })

        convert_to_micro = [
            DatadogMetrics.CONTAINER_CPU_USAGE,
            DatadogMetrics.CONTAINER_CPU_USER,
            DatadogMetrics.CONTAINER_CPU_SYSTEM,
            DatadogMetrics.CONTAINER_CPU_THROTTLED,
        ]

        for key in convert_to_micro:
            self.stats[key] *= MICRO_TO_NANOSECONDS

    def _read_memory_stats(self):

        self.stats[DatadogMetrics.CONTAINER_MEMORY_SWAP] = self._read_float(BASE_PATH, 'memory.swap.current')
        self.stats[DatadogMetrics.CONTAINER_MEMORY_USAGE] = self._read_float(BASE_PATH, 'memory.current')

        max_memory = self._read(BASE_PATH, "memory.max").rstrip()
        # https://github.com/DataDog/datadog-agent/blob/f64dc33995a71e9a4cf9e0e523b28e9466e41e82/pkg/util/cgroups/file.go#L94
        # Return None if memory.max == "max"
        self.stats[DatadogMetrics.CONTAINER_MEMORY_LIMIT] = None if max_memory == "max" else max_memory


        self._read_multi_stats(BASE_PATH, 'memory.stat', {
            DatadogMetrics.CONTAINER_MEMORY_CACHE : "file",
            DatadogMetrics.CONTAINER_MEMORY_RSS: "anon",
            "kernel_stack": "kernel_stack",
            "slab": "slab",
        })

        self.stats[DatadogMetrics.CONTAINER_MEMORY_KERNEL] = self.stats["kernel_stack"] = self.stats["slab"]

        self._read_multi_stats(BASE_PATH, 'memory.events', {
            DatadogMetrics.CONTAINER_MEMORY_OOM_EVENTS: "oom",
        })

    def collect(self):
        return self.stats
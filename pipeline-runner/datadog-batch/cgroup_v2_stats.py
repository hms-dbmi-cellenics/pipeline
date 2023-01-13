# These are compiled from:
# https://github.com/DataDog/datadog-agent/blob/main/pkg/collector/corechecks/containers/generic/processor.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv2_cpu.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv2_memory.go
# https://github.com/DataDog/datadog-agent/blob/main/pkg/util/cgroups/cgroupv2_io.go

from constants import DatadogMetrics, MetricsGroup

BASE_PATH = '/sys/fs/cgroup'
MICRO_TO_NANOSECONDS = 1000

class CGroupV2Stats:

    def __init__(self):
        self.stats = {
            MetricsGroup.CPU: {},
            MetricsGroup.MEMORY: {},
            MetricsGroup.IO: {}
        }
        self._read_cpu_stats()
        self._read_memory_stats()
        self._read_io_stats()

    def _read(self, folder, file):
        with open(f"{folder}/{file}") as f:
            return f.readline()

    def _read_float(self, folder, file):
        return float(self._read(folder, file))

    def _read_multi_stats(self, folder, file, group, mapping):
        lines = []
        value_mapping = {}

        with open(f"{folder}/{file}") as f:
            lines = f.readlines()

        for line in lines:
            key, value = line.split(" ")
            value_mapping[key] = value

        for stats_key, value_key in mapping.items():
            self.stats[group][stats_key] = float(value_mapping[value_key])

    def _read_cpu_stats(self):
        self._read_multi_stats(BASE_PATH, "cpu.stat", MetricsGroup.CPU,
        {
            DatadogMetrics.CONTAINER_CPU_USAGE: 'usage_usec',
            DatadogMetrics.CONTAINER_CPU_USER: 'user_usec',
            DatadogMetrics.CONTAINER_CPU_SYSTEM: 'system_usec',
            DatadogMetrics.CONTAINER_CPU_THROTTLED_PERIODS: 'nr_throttled',
            DatadogMetrics.CONTAINER_CPU_THROTTLED: 'throttled_usec'
        })

        convert_to_microseconds = [
            DatadogMetrics.CONTAINER_CPU_USAGE,
            DatadogMetrics.CONTAINER_CPU_USER,
            DatadogMetrics.CONTAINER_CPU_SYSTEM,
            DatadogMetrics.CONTAINER_CPU_THROTTLED,
        ]

        for key in convert_to_microseconds:
            self.stats[MetricsGroup.CPU][key] *= MICRO_TO_NANOSECONDS

    def _read_memory_stats(self):

        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_SWAP] = self._read_float(BASE_PATH, 'memory.swap.current')
        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_USAGE] = self._read_float(BASE_PATH, 'memory.current')

        max_memory = self._read(BASE_PATH, "memory.max").rstrip()
        # https://github.com/DataDog/datadog-agent/blob/f64dc33995a71e9a4cf9e0e523b28e9466e41e82/pkg/util/cgroups/file.go#L94
        # Return None if memory.max == "max"
        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_LIMIT] = None if max_memory == "max" else max_memory


        self._read_multi_stats(BASE_PATH, 'memory.stat', MetricsGroup.MEMORY,
        {
            DatadogMetrics.CONTAINER_MEMORY_CACHE : "file",
            DatadogMetrics.CONTAINER_MEMORY_RSS: "anon",
            "kernel_stack": "kernel_stack",
            "slab": "slab",
        })

        self.stats[MetricsGroup.MEMORY][DatadogMetrics.CONTAINER_MEMORY_KERNEL] = self.stats[MetricsGroup.MEMORY]["kernel_stack"] + self.stats[MetricsGroup.MEMORY]["slab"]

        self._read_multi_stats(BASE_PATH, 'memory.events', MetricsGroup.MEMORY,
        {
            DatadogMetrics.CONTAINER_MEMORY_OOM_EVENTS: "oom",
        })

    def _read_io_stats(self):
        with open(f"{BASE_PATH}/io.stat") as f:

            # File has the format:
            # 254:0 rbytes=32563200 wbytes=8192 rios=879 wios=2 dbytes=0 dios=0
            lines = f.readlines()

            for line in lines:
                device_id, rbytes, wbytes, rios, wios, dbytes, dios = line.split(" ")

                self.stats[MetricsGroup.IO][device_id] = {
                    DatadogMetrics.CONTAINER_IO_READ: float(rbytes.split("=")[1]),
                    DatadogMetrics.CONTAINER_IO_WRITE: float(wbytes.split("=")[1]),
                    DatadogMetrics.CONTAINER_IO_READ_OPERATIONS: float(rios.split("=")[1]),
                    DatadogMetrics.CONTAINER_IO_WRITE_OPERATIONS: float(wios.split("=")[1])
                }

    def collect(self):
        return self.stats
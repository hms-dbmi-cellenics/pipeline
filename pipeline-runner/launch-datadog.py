import subprocess

print("Launching Datadog subprocess...")
subprocess.Popen(["python3", "./datadog-batch/monitor.py"])
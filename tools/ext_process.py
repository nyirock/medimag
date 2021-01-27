from subprocess import Popen, PIPE


def run_external_command(cmd):
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise RuntimeError(f"process {process.args} failed with code {process.returncode}:\n {stdout} {stderr}")
    else:
        # print(stdout)
        return stdout
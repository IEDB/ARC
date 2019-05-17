import subprocess
import os

def _run_cmd(cmd, input_string=''):
        """
        Run the cmd with input_string as stdin and return output.
        """
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
               stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
               env=dict(os.environ, my_env_prop='value'), shell=True)
        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
        return out

def install(quiet):
    if quiet:
        _run_cmd("sh ARC/build_pipeline/RUN_pipeline.sh")
    else:
        print(_run_cmd("sh ARC/build_pipeline/RUN_pipeline.sh"))
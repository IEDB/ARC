import subprocess
import os
from datetime import date

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

def update(archive):
    if archive:
        today = date.today()
        folder_name = today.strftime("%d_%m_%Y_HMMs")
        abs_path = os.path.join(os.path.dirname(__file__),'../data')
        if not os.path.exists("%s/HMM_archive" % abs_path):
            print(_run_cmd(("mkdir %s/HMM_archive" % abs_path)))
        print(_run_cmd(("mv %s/HMMs %s/HMM_archive/%s") % (abs_path, abs_path, folder_name)))
        print(_run_cmd(("sh %s/build_pipeline/RUN_pipeline.sh") % os.path.dirname(__file__)))
    else:
        print(_run_cmd(("sh %s/build_pipeline/RUN_pipeline.sh") % os.path.dirname(__file__)))
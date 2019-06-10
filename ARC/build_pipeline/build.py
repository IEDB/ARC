import subprocess
import os
from datetime import date

package_directory = os.path.dirname(os.path.abspath(__file__))

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
        _run_cmd(("sh %s/RUN_pipeline.sh") % package_directory)
    else:
        print(_run_cmd(("sh %s/RUN_pipeline.sh") % package_directory))

def update(archive):
    if archive:
        today = date.today()
        folder_name = today.strftime("%d_%m_%Y_HMMs")
        if not os.path.exists("%s/../data/HMM_archive" % package_directory):
            print(_run_cmd(("mkdir %s/../data/HMM_archive" % package_directory)))
        print(_run_cmd(("mv %s/../data/HMMs %s/../data/HMM_archive/%s") % (package_directory, package_directory, folder_name)))
        print(_run_cmd(("sh %s/RUN_pipeline.sh") % package_directory))
    else:
        print(_run_cmd(("sh %s/RUN_pipeline.sh") % package_directory))
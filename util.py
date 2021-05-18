
import os
import subprocess

def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def run_external_cmd(cmd, out, cwd=None):
    """
    Run an external command/tool
    """
    command = cmd.split()
    with open(out, 'w') as fo:
        process = subprocess.run(command, stdout=fo, stderr=subprocess.PIPE, cwd=cwd)
        if process.returncode != 0:
            return process.stderr.decode()
        return None


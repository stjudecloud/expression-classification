import subprocess


def sh(command, shell=True, check=True):
  """A simple wrapper around running a shell command.

    This command will raise a RuntimeError if the command fails.
    """

  return subprocess.run(command, capture_output=True, shell=shell, check=check)

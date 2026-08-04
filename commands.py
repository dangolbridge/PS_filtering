import subprocess
from collections.abc import Sequence


def execute_command(
    command: Sequence[str],
    dry_run: bool = False,
) -> None:
    """Execute an external command and raise an error if it fails."""
    command = [str(item) for item in command]

    print(subprocess.list2cmdline(command))

    if dry_run:
        return

    subprocess.run(
        command,
        check=True,
        text=True,
    )

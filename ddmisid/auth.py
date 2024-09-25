"""Authentication utilities"""

import subprocess
import getpass


def kinit(username: str) -> None:
    """Perform kinit for Kerberos authentication with the CERN domain."""
    # Ensure the username is in the format <user>@CERN.CH
    if not username.endswith("@CERN.CH"):
        username = f"{username}@CERN.CH"

    # Prompt for the password securely
    password = getpass.getpass(prompt=f"Enter password for {username}: ")

    try:
        # Use subprocess to run kinit and pass the password via stdin
        process = subprocess.Popen(
            ["kinit", username],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate(
            input=password.encode()
        )  # Send password to kinit
        if process.returncode == 0:
            print(f"Successfully authenticated as {username}.")
        else:
            print(f"Authentication failed: {stderr.decode()}")
    except Exception as e:
        print(f"Error during kinit: {e}")

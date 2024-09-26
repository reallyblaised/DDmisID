"""Authentication utilities"""

import subprocess
import getpass
from loguru import logger


def kinit(username: str) -> None:
    """Perform kinit for Kerberos authentication with the CERN domain."""
    # Ensure the username is in the format <user>@CERN.CH
    if not username.endswith("@CERN.CH"):
        username = f"{username}@CERN.CH"

    # Allow up to three authentication attempts
    attempts = 0
    while attempts < 3:
        # Prompt for the password securely
        password = getpass.getpass(
            prompt=f"Initialising Kerberos ticket. Please enter password for {username} (Attempt {attempts + 1}/3): "
        )

        try:
            # Use subprocess to run kinit and pass the password via stdin
            process = subprocess.Popen(
                ["kinit", "-r", "7d", username],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            stdout, stderr = process.communicate(
                input=password.encode()
            )  # Send password to kinit

            # Check if kinit was successful
            if process.returncode == 0:
                logger.info(
                    f"Authentication successful. Kerberos ticket initialised for user ID {username}."
                )
                return  # Exit function on success
            else:
                attempts += 1
                logger.error(f"Authentication failed: {stderr.decode().strip()}")

        except Exception as e:
            logger.error(f"Error during kinit: {e}")
            attempts += 1

    # If all attempts fail, raise an exception
    raise RuntimeError(
        f"Failed to authenticate after 3 attempts for user ID {username}."
    )

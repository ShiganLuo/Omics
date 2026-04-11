import ftplib
import os

def download_files_from_ftp(server, username, password, remote_dir, local_dir):
    """
    Download all files from a specified directory on an FTP server to a local directory.

    :param server: FTP server address
    :param username: FTP username
    :param password: FTP password
    :param remote_dir: Directory on the FTP server to download files from
    :param local_dir: Local directory to save the downloaded files
    """
    try:
        # Connect to the FTP server
        ftp = ftplib.FTP(server)
        ftp.login(user=username, passwd=password)

        # Change to the specified remote directory
        ftp.cwd(remote_dir)

        # Ensure the local directory exists
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)

        # List files in the remote directory
        files = ftp.nlst()

        for file in files:
            local_file_path = os.path.join(local_dir, file)
            with open(local_file_path, 'wb') as local_file:
                ftp.retrbinary(f'RETR {file}', local_file.write)
                print(f"Downloaded: {file}")

        # Close the FTP connection
        ftp.quit()
        print("All files downloaded successfully.")

    except ftplib.all_errors as e:
        print(f"FTP error: {e}")

if __name__ == "__main__":
    # Example usage
    ftp_server = "ftp.example.com"
    ftp_username = "your_username"
    ftp_password = "your_password"
    remote_directory = "/path/to/remote/directory"
    local_directory = "./local_downloads"

    download_files_from_ftp(ftp_server, ftp_username, ftp_password, remote_directory, local_directory)
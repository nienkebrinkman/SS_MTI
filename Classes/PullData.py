"""
Pulls data from the MQS-data server
"""
import pysftp


class PullData:
    def __init__(self, args):
        self.args = args

    def catalog(self, cat_filepath):
        """Pull latest catalog.xml file to local computer"""

        command = ""

        pass

    def inventory(self, inv_filepath):
        """Pull latest inventory.xml file to local computer"""
        pass

    def get_waveforms(self, event, inv, remote_folder):
        with pysftp.Connection(
            host=self.args.hostname, username=self.args.username, password=self.args.password
        ) as sftp:
            print("Connection succesfully established ... ")
            if sftp.exists(remote_folder):
                event.read_waveforms(inv=inv, sc3dir=remote_folder)
            else:
                raise FileNotFoundError("Remote folder: {} does not exist")

        return event




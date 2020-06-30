import subprocess
import obspy
from os.path import exists as exist
from mqs_reports.catalog import Catalog


class GetData:
    def __init__(self):
        pass

    def read_inv(self, inv_path: str) -> None:
        """
        Read inventory file
        :param inv_path: path to inventory.xml file
        """
        self.inv = obspy.read_inventory(inv_path)

    def read_cat(self, cat_path: str) -> None:
        """
        Read catalog file
        :param inv_path: path to inventory.xml file
        :param cat_path: path to catalog.xml file
        """
        self.cat = Catalog(
            fnam_quakeml=cat_path,
            type_select=["LOW_FREQUENCY", "BROADBAND"],
            quality=["A", "B", "C"],
        )

    def mnt_remote_folder(self, host_ip: str, host_usr: str, remote_folder: str, mnt_folder: str):
        """
        Mount folder including waveform data to local machine
        :param host_ip: ip address of the server that contains waveform data
        :param host_usr: username of the host
        :param remote_folder: folder that contains waveform data
        :param mnt_folder: name of folder to mount to local machine
        """
        if not exist(mnt_folder):
            mkdir_command = "mkdir {}".format(mnt_folder)
            subprocess.call(mkdir_command, shell=True)

        mount_command = "sshfs {}@{}:{} {}".format(host_usr, host_ip, remote_folder, mnt_folder)
        subprocess.call(mount_command, shell=True)
        print("Mounted to {}".format(mnt_folder))

    def unmnt_remote_folder(self, mnt_folder):
        """
        unmount folder including waveform data to local machine
        :param mnt_folder: name of folder to mount to local machine
        """
        unmount_command = "fusermount -u {}".format(mnt_folder)
        subprocess.call(unmount_command, shell=True)
        print("UN-mounted to {}".format(mnt_folder))

    def read_events_from_cat(
        self,
        event_params: dict,
        cat: obspy.Catalog,
        inv: obspy.Inventory,
        local_folder: str,
        host_name: str = None,
        user_name: str = None,
        remote_folder: str = None,
    ) -> None:
        """
        Read mars-events from catalog.xml and adding waveform data to events.
        If waveforms are on server host_name, user_name, remote_folder should be specified
        :param event_params: list of event names for inversion
        :param cat: obspy.Catalog including the updated events
        :param local_folder: path to waveform on local machine
        :param host_name: Host IP address
        :param user_name: username of the server
        :param remote_folder: path to remote folder that contains waveform data
        """

        ## Create list of events to invert for:
        self.events = []

        ## Setup remote directory /mnt/marshost
        if host_name is not None:
            self.mnt_remote_folder(host_name, user_name, remote_folder, local_folder)

        #TODO: check if local folder exists, otherwise raise error.

        for i, v in event_params.items():
            try:
                event = cat.select(name=i).events[0]

                event.read_waveforms(inv=inv, sc3dir=local_folder)

                self.events.append(event)

            except IndexError as e:
                raise IndexError("Event {} does not exist in the catalog".format(i))

        ## Unmount:
        if host_name is not None:
            self.unmnt_remote_folder(local_folder)

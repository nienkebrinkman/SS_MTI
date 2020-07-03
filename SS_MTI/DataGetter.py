import subprocess
import obspy
from os.path import join as pjoin
from os.path import exists as exist
from os.path import isdir as isdir
from mqs_reports.catalog import Catalog
from typing import Union as _Union, Tuple as _Tuple


def read_inv(inv_path: str) -> None:
    """
    Read inventory file
    :param inv_path: path to inventory.xml file
    """
    inv = obspy.read_inventory(inv_path)
    return inv


def read_cat(cat_path: str) -> None:
    """
    Read catalog file
    :param inv_path: path to inventory.xml file
    :param cat_path: path to catalog.xml file
    """
    cat = Catalog(
        fnam_quakeml=cat_path, type_select=["LOW_FREQUENCY", "BROADBAND"], quality=["A", "B", "C"],
    )
    return cat


def unmnt_remote_folder(mnt_folder) -> None:
    """
    unmount folder including waveform data to local machine
    :param mnt_folder: name of folder to mount to local machine
    """
    unmount_command = "fusermount -u {}".format(mnt_folder)
    subprocess.call(unmount_command, shell=True)
    print("UN-mounted to {}".format(mnt_folder))


def mnt_remote_folder(host_ip: str, host_usr: str, remote_folder: str, mnt_folder: str) -> None:
    """
    Mount folder including waveform data to local machine
    :param host_ip: ip address of the server that contains waveform data
    :param host_usr: username of the host
    :param remote_folder: folder that contains waveform data
    :param mnt_folder: name of folder to mount to local machine
    """
    if not exist(mnt_folder):
        print("Create {} with writing permissions using sudo in terminal".format(mnt_folder))
        exit(1)
        # mkdir_command = "mkdir {}".format(mnt_folder)
        # subprocess.call(mkdir_command, shell=True)

    try:
        mount_command = "sshfs {}@{}:{} {}".format(host_usr, host_ip, remote_folder, mnt_folder)
        subprocess.call(mount_command, shell=True)
        print("Mounted to {}".format(mnt_folder))
    except Exception as e:
        unmnt_remote_folder(mnt_folder)
        raise e


def read_waveforms_from_saved_dir(
    self, file_path: str, event: obspy.core.event.Event = None
) -> _Union[obspy.core.event.Event, _Tuple[obspy.Stream, obspy.Stream]]:
    """
    Read in the obspy waveform data in this folder        
    :param file_path: path to waveform data
    :param event: add data to event if not None
    """
    sp_data = obspy.read(pjoin(file_path, "waveforms_SP.mseed"))
    vbb_data = obspy.read(pjoin(file_path, "waveforms_VBB.mseed"))
    if event is not None:
        event.waveforms_SP = sp_data
        event.waveforms_VBB = vbb_data
        return event
    else:
        return sp_data, vbb_data


def read_events_from_cat(
    self,
    event_params: dict,
    cat: obspy.Catalog,
    inv: obspy.Inventory,
    local_folder: str,
    host_name: str = None,
    user_name: str = None,
    remote_folder: str = None,
    save_file_name: str = None,
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
    :param save_file_name: name of .xml file to save the events, only saved when not None. (do not specify the entire folder, because this .xml file will be saved in /events/)
    """
    ## Variables that are set:
    dir_exist_name = "events"
    dir_exist = False

    ## Create empty list of events to invert for:
    events = []

    ## Setup remote directory /mnt/marshost
    event_names = [key for key in event_params.keys()]
    if all([isdir(path) for path in [pjoin(dir_exist_name, x) for x in event_names]]):
        data_exist = True
        print("The miniseed data is already saved in /events/.. and used for the inversion")
    else:
        if host_name is not None:
            mnt_remote_folder(host_name, user_name, remote_folder, local_folder)
            # TOO: check if local folder exists, otherwise raise error.

    for i, v in event_params.items():
        try:
            event = cat.select(name=i).events[0]

            # Reading the event waveforms:
            if data_exist:
                filepath = pjoin(dir_exist_name, i, "waveforms")
                event = read_waveforms_from_saved_dir(file_path=filepath, event=event)
            elif event.waveforms_VBB is not None:
                pass
            else:
                event.read_waveforms(inv=inv, sc3dir=local_folder)

            events.append(event)

        except IndexError as e:
            raise IndexError("Event {} does not exist in the catalog".format(i))

    ## Unmount:
    if host_name is not None and data_exist is False:
        unmnt_remote_folder(local_folder)

    ## save the events as a new catalog
    if save_file_name is not None:
        # event_cat = obspy.Catalog(self.events)
        # event_cat.write(
        #     filename=pjoin(dir_exist_name,save_file_name.strip(".xml")+ ".xml"), format="QUAKEML"
        # )
        pass
        # TODO: IMPLEMENT SAVE FUNCTION!

    return event


# class GetData:
#     def __init__(self):
#         pass

#     def read_inv(self, inv_path: str) -> None:
#         """
#         Read inventory file
#         :param inv_path: path to inventory.xml file
#         """
#         self.inv = obspy.read_inventory(inv_path)

#     def read_cat(self, cat_path: str) -> None:
#         """
#         Read catalog file
#         :param inv_path: path to inventory.xml file
#         :param cat_path: path to catalog.xml file
#         """
#         self.cat = Catalog(
#             fnam_quakeml=cat_path,
#             type_select=["LOW_FREQUENCY", "BROADBAND"],
#             quality=["A", "B", "C"],
#         )

#     def mnt_remote_folder(
#         self, host_ip: str, host_usr: str, remote_folder: str, mnt_folder: str
#     ) -> None:
#         """
#         Mount folder including waveform data to local machine
#         :param host_ip: ip address of the server that contains waveform data
#         :param host_usr: username of the host
#         :param remote_folder: folder that contains waveform data
#         :param mnt_folder: name of folder to mount to local machine
#         """
#         if not exist(mnt_folder):
#             print("Create {} with writing permissions using sudo in terminal".format(mnt_folder))
#             exit(1)
#             # mkdir_command = "mkdir {}".format(mnt_folder)
#             # subprocess.call(mkdir_command, shell=True)

#         try:
#             mount_command = "sshfs {}@{}:{} {}".format(
#                 host_usr, host_ip, remote_folder, mnt_folder
#             )
#             subprocess.call(mount_command, shell=True)
#             print("Mounted to {}".format(mnt_folder))
#         except Exception as e:
#             self.unmnt_remote_folder(mnt_folder)
#             raise e

#     def unmnt_remote_folder(self, mnt_folder) -> None:
#         """
#         unmount folder including waveform data to local machine
#         :param mnt_folder: name of folder to mount to local machine
#         """
#         unmount_command = "fusermount -u {}".format(mnt_folder)
#         subprocess.call(unmount_command, shell=True)
#         print("UN-mounted to {}".format(mnt_folder))

#     def read_waveforms_from_saved_dir(
#         self, file_path: str, event: obspy.core.event.Event = None
#     ) -> _Union[
#         obspy.core.event.Event, _Tuple[obspy.core.stream.Stream, obspy.core.stream.Stream]
#     ]:
#         """
#         Read in the obspy waveform data in this folder
#         :param file_path: path to waveform data
#         :param event: add data to event if not None
#         """
#         sp_data = obspy.read(pjoin(file_path, "waveforms_SP.mseed"))
#         vbb_data = obspy.read(pjoin(file_path, "waveforms_VBB.mseed"))
#         if event is not None:
#             event.waveforms_SP = sp_data
#             event.waveforms_VBB = vbb_data
#             return event
#         else:
#             return sp_data, vbb_data

#     def read_events_from_cat(
#         self,
#         event_params: dict,
#         cat: obspy.Catalog,
#         inv: obspy.Inventory,
#         local_folder: str,
#         host_name: str = None,
#         user_name: str = None,
#         remote_folder: str = None,
#         save_file_name: str = None,
#     ) -> None:
#         """
#         Read mars-events from catalog.xml and adding waveform data to events.
#         If waveforms are on server host_name, user_name, remote_folder should be specified
#         :param event_params: list of event names for inversion
#         :param cat: obspy.Catalog including the updated events
#         :param local_folder: path to waveform on local machine
#         :param host_name: Host IP address
#         :param user_name: username of the server
#         :param remote_folder: path to remote folder that contains waveform data
#         :param save_file_name: name of .xml file to save the events, only saved when not None. (do not specify the entire folder, because this .xml file will be saved in /events/)
#         """
#         ## Variables that are set:
#         dir_exist_name = "events"
#         dir_exist = False

#         ## Create empty list of events to invert for:
#         self.events = []

#         ## Setup remote directory /mnt/marshost
#         event_names = [key for key in event_params.keys()]
#         if all([isdir(path) for path in [pjoin(dir_exist_name, x) for x in event_names]]):
#             data_exist = True
#             print("The miniseed data is already saved in /events/.. and used for the inversion")
#         else:
#             if host_name is not None:
#                 self.mnt_remote_folder(host_name, user_name, remote_folder, local_folder)
#                 # TODO: check if local folder exists, otherwise raise error.

#         for i, v in event_params.items():
#             try:
#                 event = cat.select(name=i).events[0]

#                 # Reading the event waveforms:
#                 if data_exist:
#                     filepath = pjoin(dir_exist_name, i, "waveforms")
#                     event = self.read_waveforms_from_saved_dir(file_path=filepath, event=event)
#                 elif event.waveforms_VBB is not None:
#                     pass
#                 else:
#                     event.read_waveforms(inv=inv, sc3dir=local_folder)

#                 self.events.append(event)

#             except IndexError as e:
#                 raise IndexError("Event {} does not exist in the catalog".format(i))

#         ## Unmount:
#         if host_name is not None and data_exist is False:
#             self.unmnt_remote_folder(local_folder)

#         ## save the events as a new catalog
#         if save_file_name is not None:
#             # event_cat = obspy.Catalog(self.events)
#             # event_cat.write(
#             #     filename=pjoin(dir_exist_name,save_file_name.strip(".xml")+ ".xml"), format="QUAKEML"
#             # )
#             pass
#             # TODO: IMPLEMENT SAVE FUNCTION!

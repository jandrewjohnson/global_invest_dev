"""Downloader class that uses TaskGraph."""
import gzip
import logging
import os
import zipfile

import ecoshard
import retrying
import taskgraph


GZIP_BUFFER_SIZE = 2**20

LOGGER = logging.getLogger(__name__)


class TaskGraphDownloader(object):
    def __init__(
            self, download_dir, taskgraph_object_or_dir, n_workers=0):
        """Construct TaskGraphDownloader object.

        Parameters:
            download_dir (str): the base directory which files will be
                downloaded into.
            taskgraph_object_or_dir (str/TaskGraph): path to the taskgraph
                workspace database used to manage the TaskGraph object. This
                directory should not be used for any other file storage.
            n_workers (int): number of processes to use to simultaneously
                download ecoshards.
        """
        try:
            os.makedirs(download_dir)
        except OSError:
            pass
        if isinstance(taskgraph_object_or_dir, taskgraph.TaskGraph):
            LOGGER.debug('got taskgraph object')
            self.task_graph = taskgraph_object_or_dir
        else:
            LOGGER.debug('no taskgraph object, creating internal one')
            self.task_graph = taskgraph.TaskGraph(
                taskgraph_object_or_dir, n_workers)
        # this will be a dictionary indexed by ecoshard key to a dict
        # containing fields:
        #   'url': the original url
        #   'local_path': path to local file/dir
        #   'download_task': the taskgraph.Task object used to fetch the
        #         ecoshard
        self.download_dir = download_dir
        self.key_to_path_task_map = {}

    def __del__(self):
        pass

    def download_ecoshard(
            self, ecoshard_url, key, decompress='none', local_path=None):
        """Download the ecoshard given by `ecoshard_url`.

        This registers `key` as a resolved ecoshard file

        Parameters:
            ecoshard_url (str): ecoshard_url to an ecoshard file.
            key (str): key used to index the ecoshard for retrieval.
            decompress (str): if 'unzip' or 'gunzip' the "unzip" or
                "gunzip" algorithms are used to decompress the file
                downloaded by the ecoshard_url. A later `get` by `key` will
                return the directory the file was unzipped into unless
                `local_path` was defined, or in the case of GZ, the root file
                that was gunzipped.
            local_path (str): if None, the local path returned by `get` will
                be the path to the downloaded ecoshard. In the case of
                "unzip" decompression the desired local path will likely not
                be the .zip, but rather a file inside. This parameter can
                be used to add the local path with respect to the unziped
                archive that is returned on a `get`.

        Returns:
            None.

        """
        if key in self.key_to_path_task_map:
            raise ValueError(
                '"%s" is already a registered ecoshard as %s' % (
                    key, self.key_to_path_task_map[key]['url']))

        if decompress == 'none':
            local_ecoshard_path = os.path.join(
                self.download_dir, os.path.basename(ecoshard_url))
            download_task = self.task_graph.add_task(
                func=ecoshard.download_url,
                args=(ecoshard_url, local_ecoshard_path),
                target_path_list=[local_ecoshard_path],
                task_name='download %s' % local_ecoshard_path)
        elif decompress == 'gunzip':
            # ecoshard should end in .gz
            if not ecoshard_url.endswith('.gz'):
                raise ValueError(
                    'request to gunzip but %s does not end in .gz '
                    'extension' % ecoshard_url)
            local_ecoshard_path = os.path.join(
                self.download_dir, os.path.splitext(
                    os.path.basename(ecoshard_url))[0])
            download_task = self.task_graph.add_task(
                func=download_and_ungzip,
                args=(ecoshard_url, local_ecoshard_path),
                target_path_list=[local_ecoshard_path],
                task_name='download %s' % local_ecoshard_path)
        elif decompress == 'unzip':
            unzip_token_path = os.path.join(
                self.download_dir, '%s.UNZIPTOKEN' % os.path.basename(
                    ecoshard_url))
            local_ecoshard_path = self.download_dir
            download_task = self.task_graph.add_task(
                func=download_and_unzip,
                args=(ecoshard_url, self.download_dir, unzip_token_path),
                target_path_list=[unzip_token_path],
                task_name='download %s' % os.path.basename(ecoshard_url))
            download_task.join()
            if local_path:
                local_ecoshard_path = os.path.join(
                    self.download_dir, local_path)

        self.key_to_path_task_map[key] = {
            'url': ecoshard_url,
            'local_path': local_ecoshard_path,
            'download_task': download_task
        }
        LOGGER.debug('just added: %s', self.key_to_path_task_map[key])

    def get_path(self, key):
        """Return the filepath to the ecoshard referenced by key.

        If the ecoshard has not finished downloading this call will block.

        Parameters:
            key (str): key used to reference ecoshard when `download_ecoshard`
                was invoked. If `key` was not used in this way the call will
                result in an exception.

        Returns:
            filepath to ecoshard referenced by `key`.

        """
        if key not in self.key_to_path_task_map:
            raise ValueError('%s not a valid key' % key)
        local_path = self.key_to_path_task_map[key]['local_path']
        self.key_to_path_task_map[key]['download_task'].join()
        if not os.path.exists(local_path):
            raise RuntimeError('%s does not exist on disk' % local_path)
        return local_path

    def join(self):
        """Joins all downloading tasks, blocks until complete."""
        for download_map in self.key_to_path_task_map.values():
            download_map['download_task'].join()


@retrying.retry(wait_exponential_multiplier=1000, wait_exponential_max=10000)
def download_and_unzip(url, target_dir, target_token_path):
    """Download `url` to `target_dir` and touch `target_token_path`.

    Parameters:
        url (str); url to zipped file.
        target_dir (str): directory to download the `url` to and to later
            decompress the zip file to.


    """
    try:
        zipfile_path = os.path.join(target_dir, os.path.basename(url))
        LOGGER.debug('downloading %s', url)
        ecoshard.download_url(url, zipfile_path)

        LOGGER.debug('unzipping %s', zipfile_path)
        with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
            zip_ref.extractall(target_dir)

        LOGGER.debug('cleaning up %s', zipfile_path)
        os.remove(zipfile_path)

        LOGGER.debug('writing token %s', target_token_path)
        with open(target_token_path, 'w') as touchfile:
            touchfile.write(f'unzipped {zipfile_path}')
        LOGGER.debug('done with download and unzip')
    except Exception:
        LOGGER.exception('download error')
        raise


def download_and_ungzip(url, target_path):
    """Download the gzip file `url` to `target_path`.

    Downloads a gzipped file and uncompresses it to `target_path`. The original
    gzipped file is removed on successful completion of this funciton.

    Parameters:
        url (str): url to a gzipped file.
        target_path (str): local location of uncompressed gzipped file.

    Returns:
        None.

    """
    gzipfile_path = os.path.join(
        os.path.dirname(target_path), os.path.basename(url))
    ecoshard.download_url(url, gzipfile_path)

    with gzip.open(gzipfile_path, 'rb') as gzip_file:
        with open(os.path.splitext(gzipfile_path)[0], 'wb') as target_file:
            while True:
                content = gzip_file.read(GZIP_BUFFER_SIZE)
                if content:
                    target_file.write(content)
                else:
                    break

    os.remove(gzipfile_path)

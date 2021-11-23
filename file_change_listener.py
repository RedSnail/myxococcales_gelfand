import os
from abc import ABC, abstractmethod
from time import time, sleep
from threading import Thread
from sys import argv
from re import compile
from bot_utils import send_album, send_text


class ChangeListener(ABC):
    def __init__(self, path):
        self.path = path
        self.last_time = os.stat(self.path).st_mtime
        self.last_size = os.stat(self.path).st_size

    def on_change_listener(self, appendix):
        pass

    @abstractmethod
    def on_start_listener(self):
        pass

    def pass_to_listener(self, end_time):
        self.on_start_listener()
        while end_time < 0 or time() < end_time:
            next_time = os.stat(self.path).st_mtime
            if self.last_time != next_time:
                self.last_time = next_time
                new_size = os.stat(self.path).st_size
                with open(self.path, "rb") as polled_file:
                    polled_file.seek(self.last_size)
                    new_bytes = polled_file.read(new_size - self.last_size)
                    self.on_change_listener(new_bytes)
                self.last_size = new_size

            sleep(60)

    def listen(self, timelimit=-1):
        poll_start = time()
        print(poll_start + timelimit if timelimit > 0 else time())
        background_proc = Thread(target=self.pass_to_listener,
                                 name=f"{self.path} change listener",
                                 args=[poll_start + timelimit if timelimit > 0 else -1])
        background_proc.start()


class TelegramTextNotificator(ChangeListener):
    def __init__(self, path, prefix):
        super(TelegramTextNotificator, self).__init__(path)
        self.prefix = prefix

    def on_start_listener(self):
        send_text(f"starting listening to {self.path}")

    def on_change_listener(self, appendix):
        text = appendix.decode(encoding="ascii")
        cap_pattern = compile("Running blast analysis: ([0-9]*[.,][0-9])%")
        nums = cap_pattern.findall(text)
        if len(nums) != 0:
            send_text(self.prefix + " " + nums[-1] + "%")


for file in argv[1:]:
    text_listener = TelegramTextNotificator(file, file)
    text_listener.listen()




import os
from telegram import Bot, InputMediaPhoto
from telegram import ParseMode
from pathlib import Path

bot_token = os.environ["GOGA_TOKEN"]
my_id = os.environ["GOGA_ID"]
bot = Bot(token=bot_token)

HOME = Path(os.environ["MYXOCOCCALES_HOME"])


def send_album(img_paths):
    photos = list(map(lambda path: InputMediaPhoto(open(path, "rb")), img_paths))
    bot.send_media_group(my_id, photos)


def send_text(text):
    bot.sendMessage(chat_id=my_id, text=f"``{text}``", parse_mode=ParseMode.MARKDOWN)


def send_image(path):
    bot.send_photo(chat_id=my_id, photo=(HOME / path).open("rb"))

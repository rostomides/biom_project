import os


SECRET_KEY = "KJKLSJADJKSKJHASKHFKHASKJFH"

FILE_UPLOAD_PATH = os.path.abspath("static/uploads")

ALLOWED_FILES = ["CSV", "TSV"]

MAX_SIZE_FILE = 10 * 1000

# max size of accepted file
MAX_CONTENT_LENGTH = 2000000


 
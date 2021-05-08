from datetime import datetime

import pytz


def get_now_timestemp_as_string() -> str:
    return datetime.now(tz=pytz.timezone("Israel")).strftime("%Y_%m_%d__%H_%M_%S")

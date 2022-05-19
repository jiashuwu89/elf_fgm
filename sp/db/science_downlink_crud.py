from elfin.common import models
from sqlalchemy.orm import Session

# TODO: Figure out how to test these functions and if it's worth doing so


def get_first_science_downlink(db: Session):
    return db.query(models.ScienceDownlink).first()

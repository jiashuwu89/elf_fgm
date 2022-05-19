from elfin.common import db, models
from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from ..db import science_downlink_crud

router = APIRouter(
    prefix="/science_downlink",
    tags=["science_downlink"],
    responses={404: {"description": "Not found"}},
)


def get_session():
    if db.SESSIONMAKER is None:
        db.connect("testing")  # TODO: When ready, change to production

    session = db.SESSIONMAKER()

    try:
        yield session  # Dependency
    finally:
        session.close()


@router.get("/first")
def read_first_science_downlink(
    db: Session = Depends(get_session),
) -> models.ScienceDownlink:
    first_science_downlink = science_downlink_crud.get_first_science_downlink(db)
    return first_science_downlink

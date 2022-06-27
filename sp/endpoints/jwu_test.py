
from fastapi import APIRouter, Depends
import random

from ..db import science_downlink_crud

router = APIRouter(
    prefix="/jwu_test",
    tags=["jwu_test"],
    responses={404: {"description": "Not found"}},
)

@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1,10)
    b = random.randint(1,100)
    return (a,b)




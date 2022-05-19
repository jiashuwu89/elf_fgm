import logging

from fastapi import APIRouter

from ..utils import math_utils

router = APIRouter(
    prefix="/general",
    tags=["general"],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
def read_main():
    return {"msg": "Hello World, from general!!!"}


@router.get("/missions")
def missions():
    return [1, 2]


@router.get("/sum_numbers")
def sum_numbers(x: int, y: int):
    logger = logging.getLogger("sp.sum_numbers")
    logger.debug(f"sum_numbers debug: {x} {y}")
    logger.info(f"sum_numbers info: {x} {y}")
    logger.warning(f"sum_numbers warning: {x} {y}")
    return math_utils.add(x, y)

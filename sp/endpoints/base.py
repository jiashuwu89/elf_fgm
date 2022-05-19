from fastapi import APIRouter

router = APIRouter(
    responses={404: {"description": "Not found"}},
)


@router.get("/")
def read_main():
    return {"msg": "Hello World"}

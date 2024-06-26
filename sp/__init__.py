from fastapi import APIRouter, FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .endpoints import base, fgm_calib, general, science_downlink

router = APIRouter()
router.include_router(base.router)
router.include_router(general.router)
router.include_router(science_downlink.router)
router.include_router(fgm_calib.router)

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=[],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.include_router(router)

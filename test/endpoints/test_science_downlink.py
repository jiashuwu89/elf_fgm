from elfin.common import models
from fastapi.testclient import TestClient

from sp import app
from sp.db import science_downlink_crud
from sp.endpoints.science_downlink import get_session


def override_get_session():
    yield None


# Don't rely on actually connecting to the database
# https://fastapi.tiangolo.com/advanced/testing-database/
app.dependency_overrides[get_session] = override_get_session

client = TestClient(app)


# https://docs.pytest.org/en/6.2.x/monkeypatch.html
def test_read_main(monkeypatch):
    def mock_post(_):
        return models.CalculatedAttitude(mission_id=1, X=111, Y=222, Z=333)

    monkeypatch.setattr(science_downlink_crud, "get_first_science_downlink", mock_post)

    response = client.get("/science_downlink/first")
    assert response.status_code == 200
    assert response.json() == {"mission_id": 1, "X": 111, "Y": 222, "Z": 333}

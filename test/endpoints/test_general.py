from fastapi.testclient import TestClient

from sp import app

client = TestClient(app)


def test_read_main():
    response = client.get("/general")
    assert response.status_code == 200
    assert response.json() == {"msg": "Hello World, from general!!!"}


def test_missions():
    response = client.get("/general/missions")
    assert response.status_code == 200
    assert response.json() == [1, 2]


def test_sum_numbers():
    test_cases = [
        ({"x": 3, "y": 4}, 7),
        ({"x": 1, "y": 0}, 1),
        ({"x": 0, "y": 1}, 1),
        ({"x": -1, "y": 1}, 0),
    ]

    for params, expected in test_cases:
        response = client.get("/general/sum_numbers", params=params)
        assert response.status_code == 200
        assert response.json() == expected

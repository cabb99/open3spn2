def pytest_addoption(parser):
    parser.addoption("--skip-platform", action="store", default=None, help="Skip tests for a specific platform")

def pytest_configure(config):
    config.addinivalue_line("markers", "skip_platform: mark test to be skipped for a specific platform")
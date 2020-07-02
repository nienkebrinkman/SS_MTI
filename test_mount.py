try:
    pass  # mount code
except Exception as e:
    unmount()
    raise e

from .genecaller import GeneCaller

debug = False

if debug:
    import debugpy

    # Allow other computers to attach to debugpy at this IP address and port.
    debugpy.listen(("10.201.2.15", 5679))

    # Pause the program until a remote debugger is attached
    debugpy.wait_for_client()


def main():
    GeneCaller().call()


if __name__ == "__main__":
    main()

from http.server import SimpleHTTPRequestHandler
from socketserver import TCPServer
import sys


def main():
    class HttpRequestHandler(SimpleHTTPRequestHandler):
        extensions_map = {
            ""     : "application/octet-stream",
            ".html": "text/html",
            ".js"  : "application/x-javascript",
            ".wasm": "application/wasm",
        }

    port = 8000
    if len(sys.argv) > 1:
        port = int(sys.argv[1])

    with TCPServer(("localhost", port), HttpRequestHandler) as httpd:
        server_ip, server_port = httpd.server_address
        print(f"Serving HTTP on {server_ip} port {server_port} (http://{server_ip}:{server_port}/) ...", file=sys.stderr)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("Keyboard interrupt received, exiting.", file=sys.stderr)


if __name__ == "__main__":
    main()

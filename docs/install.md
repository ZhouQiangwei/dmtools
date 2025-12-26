# Install

## Build from source

```bash
git clone https://github.com/ZhouQiangwei/dmtools.git
cd dmtools
make -j$(nproc)
```

The `dmtools` binary will be in the repository root.

## Docker

```bash
docker build -t dmtools:0.1.0 .

docker run --rm dmtools:0.1.0 ./dmtools --help
```

You can mount data directories with `-v /data:/data` and pass paths inside the container.

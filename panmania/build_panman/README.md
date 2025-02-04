# To build panMAN using Dockerfile

This is tutorial to build panMAN using Dockerfile.

1. Copy the `CMakeLists.txt`, `docker`, and `install` directories to new panMAN directory.

2. Build docker image

   ```bash
   cd build
   docker built -t panman:alan
   ```

3. Go inside docker container. Working directory should be the `panman` directory.

    ```bash
    docker run -it -v $(pwd):/panman panman:alan bash
    ```

4. Run the install script.

    ```bash
    ./install/install.sh
    ```

5. Add `build/lib` to LD_LIBRARY_PATH.

    ```bash
    export LD_LIBRARY_PATH=$(pwd)/build/lib:$LD_LIBRARY_PATH
    ```

Now you should be able to run panMAN both inside and outside the docker container. I recommend keep the docker container running in a separate terminal during development so it's easier to rebuild.

Re-entering the docker container would require running `./install/reinstall_capnp.sh` to re-install capnp.

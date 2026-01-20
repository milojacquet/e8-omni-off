Generate .off files for Wythoffian E8 polytopes. Run with `cargo run -r -- [ARGS]`.

The mirrors are named `A0 A1 A2 A3 B0 B1 C M`, where `A0 B0 C` are the nodes on the length 4, 2, 1 ends respectively.

```
Usage: e8-omni-off [OPTIONS] [MIRRORS]...

Arguments:
  [MIRRORS]...  Mirrors

Options:
  -s, --single-vertex  Single vertex
  -v, --vertices       Vertex orbits
      --off-size       Estimated size of .off
  -o, --off <FILE>     Write .off
  -h, --help           Print help
  -V, --version        Print version
```
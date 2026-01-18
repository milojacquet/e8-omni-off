use crate::e8::MirrorSet;
use std::io::Write;

impl MirrorSet {
    pub fn write_off(self, writer: impl Write) {}
}

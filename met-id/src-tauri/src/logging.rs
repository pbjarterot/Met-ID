use lazy_static::lazy_static;
use log::{Level, Log, Metadata, Record};
use std::io::Write;
use std::{fs::File, sync::Mutex};

pub struct FileWriter {
    file: Mutex<File>,
}

impl Log for FileWriter {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            let mut file = self.file.lock().unwrap();
            writeln!(file, "{} - {}", record.level(), record.args()).unwrap();
        }
    }

    fn flush(&self) {}
}

lazy_static! {
    pub static ref LOGGER: FileWriter = {
        let desktop = dirs::desktop_dir().expect("Failed to get desktop directory");
        let log_path = desktop.join("metid_log.log");
        FileWriter {
            file: Mutex::new(File::create(log_path).unwrap()),
        }
    };
}

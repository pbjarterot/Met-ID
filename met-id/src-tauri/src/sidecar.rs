use serde::Serialize;
use std::fmt;
use tauri::api::process::{Command, CommandEvent};
use serde_json::to_vec;
use bincode::{serialize, deserialize};
use std::io::{Write, Read};
use std::process::Stdio;


// Custom error type
#[derive(Debug, Serialize)]
pub struct CommandError {
    message: String,
}

impl fmt::Display for CommandError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

#[tauri::command]
pub fn sidecar_function(sidecar_name: String, sidecar_arguments: Vec<String>) -> Result<String, CommandError> {
    println!("Starting processing1...");
    let (mut rx, _child) = Command::new_sidecar(sidecar_name)
        .expect("failed to create `my-sidecar` binary command")
        .args(sidecar_arguments)
        .spawn()
        .expect("Failed to spawn sidecar");
    
    let (tx, rx_data) = std::sync::mpsc::channel();
    println!("{:?}, {:?}", tx, rx_data);
    tauri::async_runtime::spawn(async move {
        // read events such as stdout
        while let Some(event) = rx.recv().await {
            if let CommandEvent::Stdout(line) = event {
                println!("{:?}", &line);
                tx.send(line).unwrap();
                
            }
        }
    });
    // Here, we collect all the stdout lines. We could also just return the first one, depending on the requirement.
    let mut output: String = String::new();
    while let Ok(line) = rx_data.recv() {
        println!("{:?}", &line);
        output.push_str(&line);
    }

    Ok(output)
}

#[tauri::command]
pub fn sidecar_function2(_sidecar_name: String, sidecar_arguments: Vec<String>) -> std::io::Result<Vec<i32>>{
    // Generate a large vector of strings
    let _binary_path = match std::env::consts::OS {
        "windows" => "metabolite-x86_64-pc-windows-msvc.exe",
        "macos" => {
            match std::env::consts::ARCH {
                "x86_64" => "metabolite_for_db-x86_64-apple-darwin",
                "aarch64" => "metabolite_for_db-aarch64-apple-darwin",
                _ => "",
            }
        },
        "linux" => "metabolite_for_db-x86_64-unknown-linux-gnu",
        _ => "",
    };
    //let data: Vec<String> = (0..100_000).map(|i| format!("string_{}", i)).collect();
    let data_bytes = to_vec(&sidecar_arguments).unwrap();  // Serialize the vector to JSON

    // Create a length-prefixed message
    let length = data_bytes.len() as u32;
    let length_bytes = serialize(&length).unwrap();

    println!("env: {:?}", std::env::current_dir());

    // Spawn the child process
    let mut child = std::process::Command::new("./pyinstaller/dist/metabolite-x86_64-pc-windows-msvc.exe")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    {
        // Get a handle to the child's stdin
        let stdin = child.stdin.as_mut().unwrap();
        // Write the length and data to the child's stdin
        stdin.write_all(&length_bytes)?;
        stdin.write_all(&data_bytes)?;
    }


    // Get a handle to the child's stdout
    let stdout = child.stdout.as_mut().unwrap();
    let mut length_buffer = [0u8; 4];
    stdout.read_exact(&mut length_buffer)?;
    let output_length: u32 = deserialize(&length_buffer).unwrap();
    
    let mut buffer = vec![0u8; output_length as usize];
    stdout.read_exact(&mut buffer)?;

    // Deserialize the received JSON to a vector of integers
    let received_data: Vec<i32> = serde_json::from_slice(&buffer).unwrap();
    println!("Received data length: {}\nData: {:?}", received_data.len(), received_data);
    

    // Wait for the child process to exit
    child.wait()?;

    Ok(received_data)

}


/* 
pub fn sidecar_function2(sidecar_name: String, sidecar_arguments: Vec<String>) -> std::io::Result<String> {
    // Serialize the vector to JSON
    let data_bytes = to_vec(&sidecar_arguments).unwrap();

    // Create a length-prefixed message
    let length = data_bytes.len() as u32;
    let length_bytes = serialize(&length).unwrap();

    // Get the current directory
    println!("env: {:?}", std::env::current_dir());

    // Execute the sidecar command
    let sidecar_command = Command::new_sidecar(&sidecar_name)?;

    // Spawn the child process
    let (mut rx, mut stdin, mut stdout) = sidecar_command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    // Write the length and data to the child's stdin
    if let Some(mut stdin) = stdin {
        stdin.write_all(&length_bytes)?;
        stdin.write_all(&data_bytes)?;
    }

    // Read the output from the child's stdout
    if let Some(mut stdout) = stdout {
        let mut length_buffer = [0u8; 4];
        stdout.read_exact(&mut length_buffer)?;
        let output_length: u32 = deserialize(&length_buffer).unwrap();

        let mut buffer = vec![0u8; output_length as usize];
        stdout.read_exact(&mut buffer)?;

        // Deserialize the received JSON to a vector of integers
        let received_data: Vec<i32> = serde_json::from_slice(&buffer).unwrap();
        println!("Received data length: {}\nData: {:?}", received_data.len(), received_data);
    }

    // Handle the events from the command execution
    while let Some(event) = rx.blocking_recv() {
        match event {
            CommandEvent::Stderr(line) => eprintln!("stderr: {}", line),
            CommandEvent::Stdout(line) => println!("stdout: {}", line),
            CommandEvent::Error(e) => return Err(std::io::Error::new(std::io::ErrorKind::Other, e)),
            CommandEvent::Terminated(_) => break,
            _ => todo!(),
        }
    }

    Ok("".to_string())
}
*/
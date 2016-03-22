task hello {
  command {
    echo 'hello World!'
  }
  output {
    File response = stdout()
  }
	runtime {
		docker: "ubuntu:latest"
	}
}


workflow test {
  call hello
}

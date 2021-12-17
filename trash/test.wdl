version 1.0

workflow wf {
  input {
    String str_var
  }
  call hello {input: str_var = str_var}
  output {
      String out = hello.out
  }
}

task hello {
  input {
    String str_var
  }

  command {
    echo ${str_var} > file.txt
  }

  # runtime {
  #   docker: "broadinstitute/my_image"
  # }
  output {
   # File out = read_lines(file.txt)
    File out = "file.txt"
  }
  # output {
  #   String out = read_string(stdout())
  # }
}
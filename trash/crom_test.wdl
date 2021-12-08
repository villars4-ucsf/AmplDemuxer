workflow myWorkflow {
    input {
        String str_var
    }
    call myTask {input: str_var = str_var}
}

task myTask {
    input {
        String str_var
    }
    command {
        echo str_var
    }
    output {
        String out = read_string(stdout())
    }
}

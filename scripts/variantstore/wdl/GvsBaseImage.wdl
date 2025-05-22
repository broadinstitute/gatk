workflow GvsBaseImage {
    String variants_docker

    call BaseImage {
        input:
            variants_docker = variants_docker,
    }

    output {
    }
}

task BaseImage {
    String variants_docker

    command <<<
        cat <<EOF
        import requests
        print(requests)
        EOF > script.py

        python script.py
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
    }
}
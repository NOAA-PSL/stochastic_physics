From noaaepic/ubuntu20.04-gnu9.3-hpc-stack:v1.2

CMD ["/bin/bash"]

ENV HOME=/home/builder

COPY --chown=builder:builder . $HOME/stochastic_physics

USER builder

ENV USER=builder

WORKDIR $HOME/stochastic_physics/unit_tests

RUN /bin/sh -c $HOME/stochastic_physics/unit_tests/build_standalone.sh

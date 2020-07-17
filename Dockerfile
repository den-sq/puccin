# Build Stage
FROM rust AS build-stage

ADD . /usr/src/puccin
WORKDIR /usr/src/myapp

RUN cargo build --release

# Final Stage
FROM scratch

ARG GIT_COMMIT
ARG VERSION
LABEL REPO="https://github.com/den-sq/puccin"
LABEL GIT_COMMIT=$GIT_COMMIT
LABEL VERSION=$VERSION

WORKDIR /usr/local/bin

COPY --from=build-stage /usr/src/puccin/bin/puccin /opt/puccin/bin/
RUN chmod +x /usr/local/bin/puccin

CMD /usr/local/bin/puccin

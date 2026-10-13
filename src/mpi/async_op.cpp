#include "async_op.h"

void AsyncOp::Request::wait() {
    if (!finished) {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        finished = true;
    }
}

void AsyncOp::Request::cancel() {
    if (!finished) {
        MPI_Cancel(&request);
        finished = true;
    }
}

bool AsyncOp::Request::checkFinished() {
    if (finished) {
        return true;
    }
    int flag = 0;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
    if (flag > 0) {
        finished = true;
    }
    return finished;
}

void AsyncOp::wait() {
    if (req_handle) {
        req_handle->wait();
        req_handle = nullptr;
    }
}

void AsyncOp::cancel() {
    if (req_handle) {
        req_handle->cancel();
        req_handle = nullptr;
    }
}

bool AsyncOp::checkFinished() {
    if (!req_handle) {
        return true;
    }
    return req_handle->checkFinished();
}

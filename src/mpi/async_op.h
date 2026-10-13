#pragma once

#include <memory>
#include <vector>
#include <iostream>

#include <mpi.h>

class AsyncOp {
public:
    class Request {
    public:
        Request(MPI_Request req) :
            request(req) {
        }
        ~Request() {
            wait();
        }

        void wait();
        void cancel();
        bool checkFinished();
    private:
        MPI_Request request;
        bool finished = false;
    };

public:
    AsyncOp() {}
    AsyncOp(MPI_Request req) {
        req_handle = std::make_shared<Request>(req);
    }
    ~AsyncOp() {}

    void wait();
    void cancel();
    bool checkFinished();
private:
    std::shared_ptr<Request> req_handle;
};

class AsyncOps {
public:
    using iterator = std::vector<AsyncOp>::iterator;
    using const_iterator = std::vector<AsyncOp>::const_iterator;
public:
    AsyncOps() {}

    void add(const AsyncOp &op) {
        _ops.push_back(op);
    }

    void add(const AsyncOps &ops) {
        if (ops.size() > 0) {
            _ops.insert(_ops.end(), ops._ops.begin(), ops._ops.end());
        }
    }

    void add(MPI_Request req) {
        add(AsyncOp(req));
    }

    void wait() {
        for (auto &op: _ops) {
            op.wait();
        }
        _ops.clear();
    }

    void cancel() {
        for (auto &op: _ops) {
            op.cancel();
        }
        _ops.clear();
    }

    size_t size() const {
        return _ops.size();
    }

    iterator begin() {
        return _ops.begin();
    }

    iterator end() {
        return _ops.end();
    }

    const_iterator begin() const {
        return _ops.begin();
    }

    const_iterator end() const {
        return _ops.end();
    }

private:
    std::vector<AsyncOp> _ops;
};

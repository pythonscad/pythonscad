/**
 * Main-thread client for the notebook Python kernel worker.
 */
import {MSG, resolveWorkerResponse,} from './notebook-protocol.mjs';

export class NotebookKernel
{
  constructor(options = {})
  {
    this.workerUrl = options.workerUrl ?? './notebook-worker.mjs';
    this.onDebug = options.onDebug ?? (() => {});
    this.worker = null;
    this.nextId = 1;
    this.pending = new Map();
    this.ready = false;
    this.busy = false;
    this.queue = [];
    this.cachedBuffers = null;
  }

  _debug(msg, type) { this.onDebug(msg, type); }

  _attachWorker(worker)
  {
    this.worker = worker;
    worker.onmessage = (event) => this._onMessage(event.data);
    worker.onerror = (event) => this._onWorkerError(event);
    worker.onmessageerror = (event) => {
      this._debug('Worker message error: ' + event, 'error');
    };
  }

  _onWorkerError(event)
  {
    const message = event.message || 'Worker script error';
    this._debug('Worker error: ' + message, 'error');
    this._failAllPending(new Error(message), true);
    this.ready = false;
  }

  _failAllPending(error, crashed = false)
  {
    error.crashed = crashed;
    for (const [, entry] of this.pending) entry.reject(error);
    this.pending.clear();
    for (const job of this.queue) job.reject(error);
    this.queue = [];
    this.busy = false;
  }

  _onMessage(data)
  {
    if (!data || typeof data.id !== 'number') return;
    resolveWorkerResponse(this.pending, data);
  }

  _request(type, payload, transfer)
  {
    return new Promise((resolve, reject) => {
      const id = this.nextId++;
      this.pending.set(id, {resolve, reject});
      try {
        this.worker.postMessage({type, id, ...payload}, transfer || []);
      } catch (e) {
        this.pending.delete(id);
        reject(e);
      }
    });
  }

  _enqueue(op)
  {
    return new Promise((resolve, reject) => {
      this.queue.push({op, resolve, reject});
      this._pumpQueue();
    });
  }

  _pumpQueue()
  {
    if (this.busy || this.queue.length === 0) return;
    const job = this.queue.shift();
    this.busy = true;
    job.op().then(job.resolve, job.reject).finally(() => {
      this.busy = false;
      this._pumpQueue();
    });
  }

  async spawn(wasmBinary, dataBinary)
  {
    this.terminate();
    if (wasmBinary && dataBinary) {
      // Keep copies for respawn: transferred buffers are neutered on the main thread.
      this.cachedBuffers = {
        wasmBinary: wasmBinary.slice(0),
        dataBinary: dataBinary.slice(0),
      };
    }
    const worker = new Worker(this.workerUrl, {type: 'module'});
    this._attachWorker(worker);

    const transfer = [];
    const payload = {};
    if (wasmBinary) {
      payload.wasmBinary = wasmBinary;
      transfer.push(wasmBinary);
    }
    if (dataBinary) {
      payload.dataBinary = dataBinary;
      transfer.push(dataBinary);
    }

    const result = await this._enqueue(() => this._request(MSG.INIT, payload, transfer));
    this.ready = Boolean(result.ready);
    if (!this.ready) throw new Error('Worker failed to initialize kernel');
    return result;
  }

  async respawnAfterCrash()
  {
    if (!this.cachedBuffers) throw new Error('No cached WASM buffers for respawn');
    this._debug('Respawning worker after crash…', 'info');
    return this.spawn(this.cachedBuffers.wasmBinary, this.cachedBuffers.dataBinary);
  }

  terminate()
  {
    if (this.worker) {
      this.worker.terminate();
      this.worker = null;
    }
    this._failAllPending(new Error('Worker terminated'), false);
    this.ready = false;
  }

  evaluate(code)
  {
    if (!this.ready || !this.worker) {
      return Promise.reject(new Error('Kernel not ready'));
    }
    return this._enqueue(() => this._request(MSG.EVALUATE, {code: String(code ?? '')}));
  }

  exportModel(format)
  {
    if (!this.ready || !this.worker) {
      return Promise.reject(new Error('Kernel not ready'));
    }
    return this._enqueue(() => this._request(MSG.EXPORT, {format: String(format || 'stl')}));
  }
}

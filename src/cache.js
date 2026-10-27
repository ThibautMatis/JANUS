// Mémoire locale (IndexedDB) des références publiques déjà téléchargées depuis NCBI.
// Une accession versionnée (ex. NG_012772.3) ne change jamais : on peut la garder indéfiniment.
// Le cache est facultatif : toute erreur (navigation privée, stockage bloqué) le rend simplement inactif.

const DB = 'svib', STORE = 'references';
let dbPromise = null;

function open() {
  if (!dbPromise) {
    dbPromise = new Promise(resolve => {
      try {
        const req = indexedDB.open(DB, 1);
        req.onupgradeneeded = () => req.result.createObjectStore(STORE);
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => resolve(null);
      } catch {
        resolve(null);
      }
    });
  }
  return dbPromise;
}

function run(mode, fn) {
  return open().then(db => new Promise(resolve => {
    if (!db) return resolve(null);
    try {
      const tx = db.transaction(STORE, mode);
      const req = fn(tx.objectStore(STORE));
      tx.oncomplete = () => resolve(req.result ?? null);
      tx.onerror = tx.onabort = () => resolve(null);
    } catch {
      resolve(null);
    }
  }));
}

export const cacheGet = key => (typeof indexedDB === 'undefined' ? Promise.resolve(null) : run('readonly', s => s.get(key)));
export const cachePut = (key, value) => (typeof indexedDB === 'undefined' ? Promise.resolve(null) : run('readwrite', s => s.put(value, key)));
export const cacheKeys = () => (typeof indexedDB === 'undefined' ? Promise.resolve([]) : run('readonly', s => s.getAllKeys()).then(k => k || []));
export const cacheClear = () => (typeof indexedDB === 'undefined' ? Promise.resolve(null) : run('readwrite', s => s.clear()));

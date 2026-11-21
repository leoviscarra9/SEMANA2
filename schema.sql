-- SIMTELE-RISK · esquema relacional del historial de ejecuciones (equivale al almacén IndexedDB "ejecuciones")
CREATE TABLE ejecucion (
  id         BIGINT PRIMARY KEY,                 -- marca de tiempo (ms)
  fecha      TEXT NOT NULL,
  modulo     TEXT NOT NULL CHECK (modulo IN ('Detención','Deterioro (sist. lineales)','Deterioro (raíces)','Oscilaciones','Sensores')),
  escenario  TEXT NOT NULL CHECK (escenario IN ('NORMAL','ESTRES','MAL')),
  parametros TEXT NOT NULL                       -- JSON con los parámetros de entrada
);
CREATE TABLE metrica (
  id           INTEGER PRIMARY KEY AUTOINCREMENT,
  ejecucion_id BIGINT NOT NULL REFERENCES ejecucion(id) ON DELETE CASCADE,
  nombre       TEXT NOT NULL,                    -- p. ej. e4, kappa, eL, eS
  valor        DOUBLE PRECISION NOT NULL
);
CREATE INDEX idx_ejecucion_modulo ON ejecucion(modulo, escenario);

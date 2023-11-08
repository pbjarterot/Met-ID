// jest.config.js or jest.config.ts
module.exports = {
  extensionsToTreatAsEsm: ['.ts'],
  transform: {
    '^.+\\.jsx?$': 'babel-jest',
    '^.+\\.tsx?$': 'ts-jest',
  },
  preset: 'ts-jest',
  testEnvironment: 'node',
  // setupFilesAfterEnv: ['<rootDir>/jest.setup.ts'], // If you have setup files
  // The glob patterns Jest uses to detect test files
  testMatch: [
    '**/__tests__/**/*.ts',
    '**/?(*.)+(spec|test).ts',
  ],
  // Module file extensions for importing
  moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
};
import { count_metabolites } from "./ms1_test_functions";

// math.test.ts
function sum(a: number, b: number): number {
  return a + b;
}

describe('sum', () => {
  it('adds 1 + 2 to equal 3', () => {
    expect(sum(1, 2)).toBe(3);
  });
});


describe('count metabolites', () => {
  it('Checks number of metabolites in database', () => {
    expect(count_metabolites("HMDB (All)", "FMP-10", ["Endogenous", "Exogenous", "Unspecified"], ["Phenols", "Primary Amines"])).toBe(26202);
  });
});
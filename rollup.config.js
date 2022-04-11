import { uglify } from "rollup-plugin-uglify"

export default {
  input: "./src/index",
  plugins: [
    uglify({
      compress: {
        pure_getters: true,
        unsafe: true,
        unsafe_comps: true,
      },
    }),
  ],
  output: {
    name: "lunar",
    file: "./dist/index.min.js",
    format: "umd",
  },
}

/** @type {import('@docusaurus/types').DocusaurusConfig} */
const path = require("path");

module.exports = {
  title: 'Genomic Variant Store',
  tagline: 'Genomic Variant Store',
  url: 'https://broadinstitute.github.io',
  baseUrl: '/gvs/', // FIXME: /gvs/
  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  trailingSlash: false,
  organizationName: 'broadinstitute', // Usually your GitHub org/user name.
  projectName: 'gvs', // Usually your repo name.
  themeConfig: {
    docs: {
      sidebar: {
        hideable: true,
      },
    },
    //  sidebar.hideable: true,
    navbar: {
      title: 'Genomic Variant Store',
      // logo: {
      //   alt: 'My Site Logo',
      //   src: 'img/logo.svg',
      // },
      items: [
        {
          type: 'doc',
          docId: 'get-started',
          position: 'right',
          label: 'Documentation',
        },
        {to: '/blog', label: 'Blog', position: 'right'},
        {
          href: 'https://github.com/broadinstitute/gatk',
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
      ],
    },
    // to integrate Hotjar feedback
    hotjar: {
      siteId: '2427684',
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Get Started',
              to: '/docs/get-started',
            },
          ],
        },
        {
          title: 'Guide and Policy',
          items: [
            {
              label: 'Contribution Guide',
              href: '/docs/contribution/README',
            },
            {
              label: 'Privacy',
              href: '/privacy',
            },
          ],
        },
        {
          title: 'Resources',
          items: [
            {
              label: 'Blog',
              to: '/blog',
            },
            {
              label: 'Changelog',
              to: 'https://github.com/broadinstitute/gatk/releases',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/broadinstitute/gatk',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Copyright © Data Sciences Platform, Broad Institute.`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        gtag: {
          trackingID: 'G-NZ4BLRHQYS',
          anonymizeIP: true,
        },
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/broadinstitute/gatk/scripts/variantstore/beta_docs/website/',
          showLastUpdateAuthor: true,
          showLastUpdateTime: true,
          sidebarCollapsible: true,
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
